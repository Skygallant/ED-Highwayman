use anyhow::{bail, Context, Result};
use kiddo::{KdTree, SquaredEuclidean};
use serde_json::Value;
use std::collections::HashMap;
use std::io::{Cursor, Read, Write};

const STORE_BYTES: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/aether.pkl"));
const GENERAL_BOOST: f32 = 6.0;

type Point = [f32; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum StarKind {
    Unknown = 0,
    General = 1,
    Neutron = 2,
}

struct InMemoryStore {
    coords: Vec<Point>,
    names: Vec<String>,
    main_stars: Vec<StarKind>,
    name_to_idx: HashMap<String, u32>,
}

fn load_store_from_bytes(bytes: &[u8]) -> Result<InMemoryStore> {
    let mut rdr = Cursor::new(bytes);
    let mut buf4 = [0u8; 4];
    rdr.read_exact(&mut buf4)?;
    let _version = u32::from_le_bytes(buf4);
    rdr.read_exact(&mut buf4)?;
    let count = u32::from_le_bytes(buf4) as usize;

    let mut coords = Vec::with_capacity(count);
    let mut names = Vec::with_capacity(count);
    let mut main_stars = Vec::with_capacity(count);

    for _ in 0..count {
        let mut buf4 = [0u8; 4];
        rdr.read_exact(&mut buf4)?;
        let x = f32::from_le_bytes(buf4);
        rdr.read_exact(&mut buf4)?;
        let y = f32::from_le_bytes(buf4);
        rdr.read_exact(&mut buf4)?;
        let z = f32::from_le_bytes(buf4);
        let mut star = [0u8; 1];
        rdr.read_exact(&mut star)?;
        rdr.read_exact(&mut buf4)?;
        let len = u32::from_le_bytes(buf4) as usize;
        let mut name_bytes = vec![0u8; len];
        rdr.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)?;
        coords.push([x, y, z]);
        names.push(name);
        let kind = match star[0] {
            1 => StarKind::General,
            2 => StarKind::Neutron,
            _ => StarKind::Unknown,
        };
        main_stars.push(kind);
    }

    let mut name_to_idx = HashMap::with_capacity(names.len());
    for (idx, name) in names.iter().enumerate() {
        name_to_idx.insert(name.clone(), idx as u32);
    }

    Ok(InMemoryStore {
        coords,
        names,
        main_stars,
        name_to_idx,
    })
}

fn build_kdtree_of(store: &InMemoryStore, kind: StarKind) -> KdTree<f32, 3> {
    let mut tree = KdTree::new();
    for (idx, coord) in store.coords.iter().enumerate() {
        if store.main_stars[idx] == kind {
            tree.add(coord, idx as u64);
        }
    }
    tree
}

fn sq_dist(a: &Point, b: &Point) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

fn find_bridge_general(
    src_idx: usize,
    dst_idx: usize,
    coords: &[Point],
    general_tree: &KdTree<f32, 3>,
    base_limit: f32,
) -> Option<(usize, f32)> {
    let general_radius_sq = base_limit.powi(2);
    let src_radius_sq = (base_limit * GENERAL_BOOST).powi(2);
    let mut best: Option<(f32, usize)> = None;

    for nn in general_tree.within::<SquaredEuclidean>(&coords[dst_idx], general_radius_sq) {
        let gen_idx = nn.item as usize;
        let dist_src_sq = sq_dist(&coords[src_idx], &coords[gen_idx]);
        if dist_src_sq > src_radius_sq {
            continue;
        }
        let dist_src = dist_src_sq.sqrt();
        let dist_dst = nn.distance.sqrt();
        let leg_base = (dist_src / GENERAL_BOOST).max(dist_dst);
        if leg_base > base_limit + 1e-6 {
            continue;
        }
        if let Some((best_base, _)) = best {
            if leg_base + 1e-6 < best_base {
                best = Some((leg_base, gen_idx));
            }
        } else {
            best = Some((leg_base, gen_idx));
        }
    }
    best.map(|(base, idx)| (idx, base))
}

fn greedy_fast_route(
    start_idx: usize,
    goal_idx: usize,
    store: &InMemoryStore,
    neutron_tree: &KdTree<f32, 3>,
    general_tree: &KdTree<f32, 3>,
    base_limit: f32,
) -> Option<(Vec<usize>, f32, usize)> {
    let mut limit = base_limit;
    let mut current = start_idx;
    let mut generals: Vec<usize> = Vec::new();
    let mut required_base = 0.0_f32;
    let mut steps = 0usize;
    let mut no_progress_iters = 0usize;

    let goal_dist = |idx: usize| -> f32 {
        sq_dist(&store.coords[idx], &store.coords[goal_idx]).sqrt()
    };

    while current != goal_idx {
        let cur_goal = goal_dist(current);
        let edge_limit_sq = (limit * (GENERAL_BOOST + 1.0_f32)).powi(2);
        let mut candidates: Vec<(f32, usize)> = Vec::new();
        for nn in neutron_tree.within::<SquaredEuclidean>(&store.coords[current], edge_limit_sq) {
            let nidx = nn.item as usize;
            if nidx == current || !matches!(store.main_stars[nidx], StarKind::Neutron) {
                continue;
            }
            let dist_goal = goal_dist(nidx);
            if dist_goal + 1e-3 < cur_goal {
                candidates.push((dist_goal, nidx));
            }
        }
        candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let mut moved = false;
        for (_, cand) in candidates {
            if let Some((bridge_general, leg_base)) =
                find_bridge_general(current, cand, &store.coords, general_tree, limit)
            {
                required_base = required_base.max(leg_base);
                generals.push(bridge_general);
                current = cand;
                steps += 1;
                moved = true;
                no_progress_iters = 0;
                break;
            }
        }

        if !moved {
            limit += 1.0;
            no_progress_iters += 1;
            if no_progress_iters % 25 == 0 {
                eprintln!(
                    "Raised limit to {:.3} after {} stalls (at {})",
                    limit, no_progress_iters, store.names[current]
                );
            }
            if limit > 1000.0 {
                eprintln!("Exceeded 1000 ly limit; aborting");
                return None;
            }
        }
    }

    Some((generals, required_base, steps.saturating_mul(2)))
}

fn load_aliases(path: &str) -> Result<HashMap<String, String>> {
    if !std::path::Path::new(path).exists() {
        let default = serde_json::json!({
            "Sol": "Jackson's Lighthouse",
            "Colonia": "Magellan"
        });
        std::fs::write(path, serde_json::to_string_pretty(&default)?)
            .context("failed to write default jumppoints.json")?;
    }
    let file = std::fs::File::open(path)?;
    let val: Value = serde_json::from_reader(file)?;
    let obj = val
        .as_object()
        .context("alias file must be a JSON object of name->system")?;
    let mut out = HashMap::new();
    for (k, v) in obj {
        if let Some(s) = v.as_str() {
            out.insert(k.clone(), s.to_string());
        }
    }
    Ok(out)
}

fn resolve_alias(name: &str, aliases: &HashMap<String, String>) -> String {
    if let Some(stripped) = name.strip_prefix("JP:") {
        if let Some(target) = aliases.get(stripped) {
            return target.clone();
        }
    }
    name.to_string()
}

fn print_route(
    start_idx: usize,
    goal_idx: usize,
    generals: &[usize],
    required_base: f32,
    store: &InMemoryStore,
) {
    let mut out = String::new();
    out.push_str(&format!(
        "Minimum jump range required: {:.2}\n",
        required_base
    ));
    out.push_str(&format!("Total neutron jumps: {}\n", generals.len()));
    out.push_str(&store.names[start_idx]);
    out.push('\n');
    for g in generals {
        out.push_str(&store.names[*g]);
        out.push('\n');
    }
    out.push_str(&store.names[goal_idx]);
    std::fs::write("route.txt", out).expect("failed to write route.txt");
}

fn main() -> Result<()> {
    println!("Reminder: only enter neutron stars for start and destination.");
    println!("Reminder: invoke custom neutron star names with the prefix JP:");

    let mut input = String::new();
    print!("Start system: ");
    std::io::stdout().flush().ok();
    std::io::stdin().read_line(&mut input)?;
    let start_name_raw = input.trim().to_string();
    input.clear();

    print!("Destination system: ");
    std::io::stdout().flush().ok();
    std::io::stdin().read_line(&mut input)?;
    let goal_name_raw = input.trim().to_string();
    input.clear();

    print!("Max jump distance (ly) [default 30]: ");
    std::io::stdout().flush().ok();
    std::io::stdin().read_line(&mut input)?;
    let mut max_distance: f32 = 30.0;
    if let Ok(v) = input.trim().parse() {
        max_distance = v;
    }

    let aliases = load_aliases("jumppoints.json").unwrap_or_default();
    let start_name = resolve_alias(&start_name_raw, &aliases);
    let goal_name = resolve_alias(&goal_name_raw, &aliases);

    let store = load_store_from_bytes(STORE_BYTES).context("failed to load embedded store")?;
    let neutron_tree = build_kdtree_of(&store, StarKind::Neutron);
    let general_tree = build_kdtree_of(&store, StarKind::General);

    let start_idx = *store
        .name_to_idx
        .get(&start_name)
        .ok_or_else(|| anyhow::anyhow!("Start system not found: {}", start_name))? as usize;
    let goal_idx = *store
        .name_to_idx
        .get(&goal_name)
        .ok_or_else(|| anyhow::anyhow!("Destination system not found: {}", goal_name))? as usize;

    if !matches!(store.main_stars.get(start_idx), Some(StarKind::Neutron)) {
        bail!("Start system must be a neutron star");
    }
    if !matches!(store.main_stars.get(goal_idx), Some(StarKind::Neutron)) {
        bail!("Destination system must be a neutron star");
    }
    if start_idx == goal_idx {
        println!("{}", store.names[start_idx]);
        println!("{}", store.names[goal_idx]);
        println!("Total general stars: 0");
        return Ok(());
    }

    eprintln!(
        "Routing {} -> {} starting at base {:.2} ly (boosted {:.2} ly)...",
        start_name,
        goal_name,
        max_distance,
        max_distance * GENERAL_BOOST
    );

    if let Some((generals, required_base, jumps)) = greedy_fast_route(
        start_idx,
        goal_idx,
        &store,
        &neutron_tree,
        &general_tree,
        max_distance,
    ) {
        eprintln!(
            "Route found with {} legs, max base {:.3} ly, jumps {}",
            generals.len(),
            required_base,
            jumps
        );
        print_route(start_idx, goal_idx, &generals, required_base, &store);
    } else {
        bail!("No route found within search bounds");
    }

    print!("Press Enter to exit...");
    std::io::stdout().flush().ok();
    let mut dummy = String::new();
    std::io::stdin().read_line(&mut dummy).ok();

    Ok(())
}
