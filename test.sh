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
../roda.exe test1.inp > out_files/output1.log
# Compare the output with expected results
# You can use `diff
echo "################Test-1####################"
echo "Integral of he2 from wfx (single center)"
diff he2.out expected_he2.out 

#------------------test-2------------------------------#
#-----------Total integral of H2O----------------------#      
../roda.exe test2.inp > out_files/output2.log
# Compare the output with expected results
# You can use `diff  
echo "################Test-2####################"
echo "Integral of water from wfx (single center)"   
diff water.out tests/expected_water.out 

#------------------test-3------------------------------#
#-----------Radial plot of He2-------------------------# 
../roda.exe test3.inp > out_files/output3.log
# Compare the output with expected results
# You can use `diff
echo "################Test-3####################"  
echo "radial plot of he2 from wfx" 
diff he2.rad expected_he2.rad >> out

#------------------test-4------------------------------#
#-----------Radial plot of H2O-------------------------# 
../roda.exe test4.inp > out_files/output4.log
# Compare the output with expected results
# You can use `diff      
echo "################Test-4####################"       
echo "radial plot of water from wfx"                         
diff water.rad expected_water.rad >> out

#------------------test-5------------------------------#
#-----------Total integral of H3 (from fchk)-----------# 
../roda.exe test5.inp > out_files/output5.log
# Compare the output with expected results
# You can use `diff                               
echo "################Test-5####################"       
echo "Integral of H3 from fchk (single center)"
diff h3.out expected_h3.out >> out

#------------------test-6------------------------------#
#-----------Total integral of Li (from fchk)-----------# 
../roda.exe test6.inp > out_files/output6.log
# Compare the output with expected results
# You can use `diff
echo "################Test-6####################"       
echo "Integral of Li from fchk"             
diff Li.out expected_Li.out 

#------------------test-7------------------------------#
#--------Multicenter integral of He2(manual quadrature)#
../roda.exe test7.inp > out_files/output7.log
# Compare the output with expected results
echo "################Test-7####################"       
echo "Multicenter integral of He2 from wfx (manual quadrature)"
diff he2_multi.out expected_he2_multi.out

#------------------test-8------------------------------#
#--------Multicenter integral of He2(automatic)--------#
../roda.exe test8.inp > out_files/output8.log
# Compare the output with expected results
echo "################Test-8####################"       
echo "Multicenter integral of He2 from wfx (automatic quadrature)"
diff he2_auto.out expected_he2_auto.out

#------------------test-9------------------------------#
#--------Multicenter integral of H2O(automatic)--------#
../roda.exe test9.inp > out_files/output9.log
# Compare the output with expected results
echo "################Test-9####################"
echo "Multicenter integral of water from wfx (automatic quadrature)"
diff water_auto.out expected_water_auto.out

#------------------test-10------------------------------#
#--------Intracule at zero of Li------------------------#
../roda.exe test10.inp > out_files/output10.log
# Compare the output with expected results
echo "################Test-10####################"
echo "Intracule at zero of Li from fchk"
diff Li_intracule_zero.out expected_Li_intracule_zero.out

#----------------Add open shell molecules!!!!----------#

echo "All tests passed."

