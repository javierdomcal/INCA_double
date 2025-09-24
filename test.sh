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
echo "Integral of he2 from wfx (single center)"
diff he2.out expected_he2.out > out

#------------------test-2------------------------------#
#-----------Total integral of H2O----------------------#      
../roda.exe test2.inp > output2.log
# Compare the output with expected results
# You can use `diff  
echo "################Test-2####################"
echo "Integral of water from wfx (single center)"   
diff water.out tests/expected_water.out >> out

#------------------test-3------------------------------#
#-----------Radial plot of He2-------------------------# 
../roda.exe test3.inp > output3.log
# Compare the output with expected results
# You can use `diff
echo "################Test-3####################"  
echo "radial plot of he2 from wfx" 
diff he2.rad expected_he2.rad >> out

#------------------test-4------------------------------#
#-----------Radial plot of H2O-------------------------# 
../roda.exe test4.inp > output4.log
# Compare the output with expected results
# You can use `diff      
echo "################Test-4####################"       
echo "radial plot of water from wfx"                         
diff water.rad expected_water.rad >> out

#------------------test-5------------------------------#
#-----------Total integral of H3 (from fchk)-----------# 
../roda.exe test5.inp > output5.log
# Compare the output with expected results
# You can use `diff                               
echo "################Test-5####################"       
echo "Integral of H3 from fchk (single center)"
diff h3.out expected_h3.out >> out

#------------------test-6------------------------------#
#-----------Total integral of Li (from fchk)-----------# 
../roda.exe test6.inp > output6.log
# Compare the output with expected results
# You can use `diff
echo "################Test-6####################"       
echo "Integral of Li from fchk"             
diff Li.out expected_Li.out >> out


#----------------Add open shell molecules!!!!----------#
cat out

echo "All tests passed."

