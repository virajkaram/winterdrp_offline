# winterdrp_offline
Some scripts to reduce winter data.

First, grab some unpacked and split images from the winter machine at caltech - `/data/loki/data/winter/<night>/raw_unpacked`.
Then, grab an observing log from the winter machine - `/data/loki/observing_logs/observation_log_<night>.csv`
In `winterdrp_offline/run.py`, edit the USER-SPECIFIABLE options at the top of the file. 
Run using `python winterdrp_offline/run.py`.
