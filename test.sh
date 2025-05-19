#!/bin/bash
#set -e  # Exit on error
#set -x  # Debug print

# Compile the program
make clean
make

# Run the program (adjust the input file name as needed)
# go to the test directory to execute the program
cd tests/
#------------------test-1------------------------------#
#-----------Total integral of He2----------------------#
../roda.exe test1.inp > output1.log
# Compare the output with expected results
# You can use `diff
echo "################Test-1####################"
diff he2.out expected_he2.out > out

#------------------test-2------------------------------#
#-----------Total integral of H2O----------------------# 
echo "################Test-2####################"       
../roda.exe test2.inp > output2.log
# Compare the output with expected results
# You can use `diff   
diff water.out tests/expected_water.out >> out

#------------------test-3------------------------------#
#-----------Radial plot of He2-------------------------#
echo "################Test-3####################"       
../roda.exe test3.inp > output3.log
# Compare the output with expected results
# You can use `diff                               
diff he2.rad expected_he2.rad >> out

#------------------test-4------------------------------#
#-----------Radial plot of H2O-------------------------# 
echo "################Test-4####################"       
../roda.exe test4.inp > output4.log
# Compare the output with expected results
# You can use `diff                               
diff water.rad expected_water.rad >> out

#----------------Add open shell molecules!!!!----------#
cat out

echo "All tests passed."

