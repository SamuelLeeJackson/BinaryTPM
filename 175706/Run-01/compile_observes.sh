START=0
END=400
STEP=10
 
for (( c=$START; c<=$END; c=c+STEP ))
do
    g++ observe_generator.cpp -DT_INERTIA=$c -std=c++1y -o observe_${c} -O3
done
