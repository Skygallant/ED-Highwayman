# The Highwayman Neutron Route Plotter for Elite: Dangerous

Download the .exe from the releases section and double click it to start plotting! ✨\
**Please note that the program needs about 5GB of free RAM to run properly.**

### Features:
- Currently the code is setup for the Caspian only.
  - (this may change in the future)
- Assign custom names to your favourite neutron stars with the jumppoints.json
  - The format is: `"CUSTOM NAME": "IN-GAME NAME"`
  - Use your custom names in the program like this: `JP:CUSTOM NAME` (case sensitive)
- The program outputs the route to a handy dandy route.txt for your perusal.

### How to use:
- Set your galmap to neutron star mode, then punch in the next system in the list as a destination.
  - The galmap should detect a nuetron star to use to get to the star you set as your destination.
  - (All the systems in the route.txt are fuel stars)

### How it works:
- There is an aether.pkl file in the repo which is a custom database for this program parsed from galaxy.json.gz from spansh.co.uk (thank you Spansh!)
- When you run the program, it searchs for a series of neutron -> fuel star -> neutron hops so that the star types alternate for simplicity.
- The galmap's neutron routing feature is utilisied so that the route.txt doesn't need to store all the jumps, only the fuel stars of each neutron -> fuel star hop.

I hope you enjoy this program and that it is of some help to people! :D