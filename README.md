# run-event-classifier
Scripts to run Akash's event classifier on PNSN events in real-time.

run_Automated_Surface_Event_Detection.py is based on this notebook: https://github.com/Akashkharita/Surface_Event_Detection/blob/main/Notebooks/Automated_Surface_Event_Detection.ipynb
I modified it so that the only input argument is an evid.  
It uses the evid to make queries to find the 10 closest stations with picks.  This works for Binder events, Jiggled events and subnet triggers.
It runs the model, then generates a figure and summarizes station and event-level results.

Results are currently being dumped to: https://seismo.ess.washington.edu/~ahutko/EVENT_CLASSIFICATION/

Currently I am only running model P_50_100_F_1_10_50 from July 3.  Akash has since updated this.

get_and_check_evids.py is the script I use to hunt down evids to feed into a cronscript.
