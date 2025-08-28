%% IEEE-5 bus system
nbus = 5;
Sbase = 100; %% MVA, Un=230kV
wb = 2*pi*50; %% rad/s
START_BUS = 0;

%%      |  From_Bus |  T_ Bus  |   R_pu    |   X_pu     |     B_pu
branch_data=[
                0         1      0.00281    0.02810    0.00593
                0         3      0.00304    0.03040    0.00548
                0         4      0.00064    0.00640    0.02605
                1         2      0.00108    0.01080    0.01543
                2         3      0.00297    0.02970    0.00562
                3         4      0.00297    0.02970    0.00562
             ];

%%         |  at_Bus   |Bus_type  | P_Load_pu | Q_Load_pu  |LF_Voltage_pu|Nominal_U_kv |
bus_data=[
                1            1     3.00000     0.98610    1.0000      230.0
                2            2     3.00000     0.98610    1.0000      230.0
                3            3     4.00000     1.31470    1.0000      230.0
             ];

      %buses with 0 load demand are not listed here
%%           |  at_Bus   | Bus_type  | Voltage_pu    |P_gen_pu   | Q_gen_pu   |Sub-transient|
gen_data=[
                3     3       1.0000     0.000000   0.00000     0.25
                0     2       1.0000     0.400000   0.00000     0.25
                2     2       1.0000     3.234900   0.00000     0.25
                4     2       1.0000     4.665100   0.00000     0.25
                0     2       1.0000     1.700000   0.00000     0.25
             ];
