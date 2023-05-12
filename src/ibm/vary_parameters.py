#!/usr/bin/env python3

import datetime as dt

exe_name ="./ipd_simulation.exe"

max_time = 50000

init_x = [0.5]
init_y = [0.5]

mu_y = [ 0.02 ]

mu_xp = 0.02
mu_yp = 0.02

nrep = 10

counter = 0

basename = "data_ipd_"
current_date = dt.datetime.now()
basename += current_date.strftime("%Y%m%d_%H%M%S%f")

for rep_i in range(0,nrep):
    for init_x_i in init_x:
        for init_y_i in init_y:
            for mu_y_i in mu_y:
                print(f"{exe_name}" +\
                        f" {max_time} {init_x_i}" +\
                        f" {init_y_i} {mu_y_i} {mu_xp} {mu_yp} {basename}_{counter}")

                counter += 1
