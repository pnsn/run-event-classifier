# run-event-classifier
(Unpolished) scripts to run Akash's event classifier on PNSN events in real-time.

run_Automated_Surface_Event_Detection.py is based on this notebook: https://github.com/Akashkharita/Surface_Event_Detection/blob/main/Notebooks/Automated_Surface_Event_Detection.ipynb
I modified it so that the only input argument is an evid.  
It uses the evid to make queries to find the 10 closest stations with picks.  This works for Binder events, Jiggled events and subnet triggers.
It runs the model, then generates a figure and summarizes station and event-level results.

Results are currently being dumped to: https://seismo.ess.washington.edu/~ahutko/EVENT_CLASSIFICATION/

Currently I am only running model P_50_100_F_1_10_50 from July 3.  Akash has since updated this.

get_and_check_evids.py is the script I use to hunt down evids to feed into a cronscript.


<b>Caveats:</b>

-Akash looked into this and his model works best for shallow events, i.e. it can get deeper events (>30km) wrong and call them EX up to 20% of the time.  Same is true for larger events M>3 where the larger the event, the more likely it is to get it wrong- likely due to longer coda which make it look like a surface event.

-The model was trained on vertical broadband and some short period data.  I apply it to only vertical channels, but I include strong-motion-only stations if they part of the 10 closest/earliest picks.

-This is V0beta of it all, feedback welcome, but right now it's just "something" for me so that I can get a feel for how this performs in the wild.  It's not meant to be a polished product.

-The model and base code that's in the repo is constantly evolving, so it looks quite messy.

<b>Reading a figure</b>:

-Bottom panel shows the average of the 10 (or however many) stations analyzed.

-The analysis window is 150 sec long, the dots located in the middle of that window.  The color of the dot corresponds to the event type with the maximum probability.

-The data plotted are 1-10Hz filtered which is the preprocessing that goes into the ML model.

-The solid black vertical line is the origin time or earliset pick time for subnet triggers.

-Above each seismogram the maximum probabilities for each event type (for any time window) is given.  "Prediction" is the maximum of these unless the non-noise probs are low.

