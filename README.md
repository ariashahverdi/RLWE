## In the following I am explaining how to run the code

### The requirements are *sage* and *protobuf*

The main code which runs the attack is **leakage_attack.py**

* You can run the challenges by using the following command

`
./run_all.sh challenges/chall-id0001-rlwed-m256-q769-l3-short-toy '1-7mod16' 'Challenges' 'Gaussian'
`

* You can run the attack on the NewHope parameter set you can use the following command for Gaussian or Binomial case

For Gaussian Distribution
`
./run_all.sh '' '1mod8' 'NewHope' 'Gaussian'
`

For Binomial Distribution
`
./run_all.sh '' '1mod8' 'NewHope' 'Binomial'
`
