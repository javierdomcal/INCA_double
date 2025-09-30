#!/bin/bash
set -e

# Compile program
make clean && make

cd tests/

# Define systems
systems=("He2" "H2O" "H3" "Li")

for sys in "${systems[@]}"; do
    echo "============================================"
    echo " Running tests for $sys"
    echo "============================================"

    cd "$sys"
    mkdir -p generated/logs generated/results
    # Run all inputs for this system
    for inp in inputs/*.inp; do
        testname=$(basename "$inp" .inp)
        echo ">>> Running $sys : $testname"
        mv aux/* .
        ../../roda.exe "$inp" > "generated/logs/${testname}.log"
        mv *.dm2p aux/
        if ls *.fchk 1> /dev/null 2>&1; then mv *.fchk aux/; fi
        if ls *.wfx 1> /dev/null 2>&1; then mv *.wfx aux/; fi

        resultfile=$(grep '\.out$' "$inp" | tail -n 1 | xargs)
        if [ -z "$resultfile" ]; then 
            echo "⚠️ Could not detect result file in $inp" 
        else 
            # Move result file to results/ folder 
            if [ -f "$resultfile" ]; then 
                mv "$resultfile" "generated/results/$resultfile" 
            fi
            # Compare each output with reference
            if [ -f "reference/$resultfile" ]; then
                if grep -q "Radial_integral error=" "generated/results/$resultfile"; then
                    echo "Checking $resultfile against expected_$resultfile"
                    gen_err=$(grep "Radial_integral error=" "generated/results/$resultfile") 
                    ref_err=$(grep "Radial_integral error=" "reference/$resultfile")
                    if [ "$gen_err" == "$ref_err" ]; then
                        echo "✅ Radial integral error matches"
                    else
                        echo "❌ Radial integral error mismatch"
                        echo "   Generated: $gen_err"
                        echo "   Expected:  $ref_err"
                    fi
                else
                    echo "No radial integral error found in $resultfile"        
                    diff "generated/results/$resultfile" "reference/$resultfile" \
                    || echo "❌ Mismatch in $sys/$testname ($resultfile)"
                fi    
            else
                echo "⚠️ No reference file for $resultfile"
            fi
        fi    
    done
    cd ..
done

echo "✅ All tests completed."
