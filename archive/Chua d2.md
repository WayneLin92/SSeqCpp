# (125, 15)
## Chua -> weinan_res
x_(125,15,0) = [1, 1, 0, 0, 1]
x_(125,15,1) = [0, 0, 0, 0, 1]
x_(125,15,2) = [0, 1, 1, 0, 0]
x_(125,15,3) = [1, 1, 1, 0, 1]
x_(125,15,4) = [1, 0, 0, 1, 1]

## Chua -> weinan_res
x_(125,15,0) = [0,1,4]
x_(125,15,1) = [4]
x_(125,15,2) = [1,2]
x_(125,15,3) = [0,1,2,4]
x_(125,15,4) = [0,3,4]

## weinan_basis -> weinan_res
(125, 15)
0, 1
1, 0,1
2, 2
3, 3
4, 4

## weinan_res -> weinan_basis
[0] -> [0,1]
[1] -> [0]
[2] -> [2]
[3] -> [3]
[4] -> [4]

# (124, 17)
## Chua -> weinan_res
x_(124,17,0) = [0, 0, 1, 1]
x_(124,17,1) = [0, 0, 0, 1]
x_(124,17,2) = [1, 1, 1, 1]
x_(124,17,3) = [0, 1, 0, 0]

## Chua -> weinan_res
x_(124,17,0) = [2,3]
x_(124,17,1) = [3]
x_(124,17,2) = [0,1,2,3]
x_(124,17,3) = [1]

## weinan_basis -> weinan_res
0, 0
1, 1
2, 3
3, 2

## weinan_res -> weinan_basis
[0] -> [0]
[1] -> [1]
[2] -> [3]
[3] -> [2]

# d2
## Chua d2
[0, 0, 0, 0]
[0, 0, 0, 0]
[1, 1, 0, 0]
[0, 0, 0, 0]
[1, 1, 0, 1]

## Chua d2
d2[0] = []
d2[1] = []
d2[2] = [0,1]
d2[3] = []
d2[4] = [0,1,3]

## weinan_res d2
[0] -> [0,1,4] -> [1,4] => []
[1] -> [4] -> [4] => []
[2] -> [1,2] -> [0,2] => [0,1] -> [2] -> [3]
[3] -> [0,1,2,4] -> [1,2,4] => []
[4] -> [0,3,4] -> [0,1,3,4] => [0,1,3] -> [1,2] -> [1,3]

[1,4] => []
[4] => []
[0,2] => [3]
[1,2,4] => []
[0,1,3,4] => [1,3]

[1] => []
[4] => []
[0] => [3]
[2] => []
[0,3] => [1,3]

[0] => [3]
[1] => []
[2] => []
[3] => [1]
[4] => []

# (124, 15)
## Chua -> weinan_res
x_(124,15,0) = [1, 1, 1, 1]
x_(124,15,1) = [1, 1, 0, 1]
x_(124,15,2) = [0, 1, 0, 0]
x_(124,15,3) = [1, 1, 0, 0]

## Chua -> weinan_res
x_(124,15,0) = [0,1,2,3]
x_(124,15,1) = [0,1,3]
x_(124,15,2) = [1]
x_(124,15,3) = [0,1]

## weinan_basis -> weinan_res
0, 1
1, 0
2, 3
3, 2

## weinan_res -> weinan_basis
[0] -> [1]
[1] -> [0]
[2] -> [3]
[3] -> [2]

# (123, 17)
## Chua -> weinan_res
x_(123,17,0) = [1, 0, 1, 1]
x_(123,17,1) = [0, 1, 1, 1]
x_(123,17,2) = [0, 0, 1, 1]
x_(123,17,3) = [1, 1, 1, 0]

## Chua -> weinan_res
x_(123,17,0) = [0,2,3]
x_(123,17,1) = [1,2,3]
x_(123,17,2) = [2,3]
x_(123,17,3) = [0,1,2]

## weinan_basis -> weinan_res
0, 1
1, 0
2, 2
3, 3

## weinan_res -> weinan_basis
[0] -> [1]
[1] -> [0]
[2] -> [2]
[3] -> [3]

# d2
## Chua d2
0,[1, 1, 0, 1]
1,[1, 1, 1, 1]
2,[0, 0, 0, 0]
3,[0, 0, 1, 0]

## Chua d2
d2[0] = [0,1,3]
d2[1] = [0,1,2,3]
d2[2] = []
d2[3] = [2]

## weinan_res d2
[0] -> [0,1,2,3] -> [0,1,2,3] => [0,1,3] -> [2] -> [2]
[1] -> [0,1,3] -> [0,1,2] => [0,1,2,3] -> [3] -> [3]
[2] -> [1] -> [0] => []
[3] -> [0,1] -> [0,1] => [2] -> [2,3] -> [2, 3]

[0,1,2,3] => [2]
[0,1,2] => [3]
[0] => []
[0,1] => [2,3]

[3] => [2,3]
[2] => [2]
[0] => []
[1] => [2,3]

[0] => []
[1] => [2,3]
[2] => [2]
[3] => [2,3]

###############
x_(101,11,0) = [1, 0, 0]
x_(101,11,1) = [0, 1, 0]
x_(101,11,2) = [0, 1, 1]

x_(101,11,0) = 0 -> 0,1
x_(101,11,1) = 1 -> 0
x_(101,11,2) = 1,2 -> 0,2

0 -> 1
1 -> 0,1
2 -> 2

0 -> 0,1
1 -> 0
2 -> 2