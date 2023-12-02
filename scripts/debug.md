deg_x=(94, 8) deg_dx=(93, 11) deg_y=(70, 6) deg_xy=(164, 14) deg_dxy=(163, 17)
x=x_{190} dx=x_1x_{189} y=v_{70} dy=0 xy=x_{868}v_0 dxy=x_1x_3x_{764}v_0
x=[0] dx=[0] y=[0] dy=[-1] xy=[6] dxy=[2]

Debug
```
del debug && xcopy auto-2023-10-22-00-36 debug && del debug\deduction.log && ss add_diff S0 164 14 3 0 0 debug && ss add_diff S0 38 4 3 0 0 debug && ss add_diff C2h6 155 9 3 0 "" debug
ss add_diff CW_2_eta 70 6 3 0 "" debug && ss add_diff S0 94 8 3 0 0 debug

ss add_diff S0 38 4 3 0 0 debug
ss add_diff S0 94 8 3 0 0 debug
 
ss add_diff C2h6 155 9 3 0 "" debug
ss add_diff C2h6 126 10 3 4 5 debug

ss add_diff C2h6 126 10 3 4 "" debug

```

```
add_diff_from_log 1087342
```

```
manual - C2h6 (126, 10) d_3[4]=[]
Error(0x75989376) - No source for the image. C2h6 deg_dx=(163, 17), dx=[4], r=2
```

Joker (163,13) d_3[1]=[6]

CW_2_eta (70, 6) d3(h0Q3[3]) = 0
S0 (94, 8) d3(x_{94, 8}) = h1x_{92, 10}
x_{94, 8} * h0Q3[3] = x_{164, 14, 2}[0]
h1x_{92, 10} * h0Q3[3] = h1h3x_{155,15}[0]
因此 CW_2_eta d3(x_{164, 14, 2}[0]) = h1h3x_{155,15}[0]

e1 * Delta^2h3^2[64] + Ph1 * d0^2h6[64]