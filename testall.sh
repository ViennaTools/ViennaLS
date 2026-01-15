rm -rf build
cmake -B build -DVIENNALS_BUILD_TESTS=ON -DVIENNALS_BUILD_EXAMPLES=ON
cmake --build build -j 16
ctest -E "Benchmark|Performance" --test-dir build
cd build/examples/AirGapDeposition
echo "=== AirGapDeposition ==="
./AirGapDeposition
cd ../Deposition 
echo "=== Deposition ==="
./Deposition
cd ../Epitaxy
echo "=== Epitaxy ==="
./Epitaxy 
cd ../GeometricAdvection
echo "=== GeometricAdvection ==="
./GeometricAdvection 
cd ../PatternedSubstrate
echo "=== PatternedSubstrate ==="
./PatternedSubstrate 
cd ../PeriodicBoundary
echo "=== PeriodicBoundary ==="
./PeriodicBoundary 
cd ../SharedLib
echo "=== SharedLib ==="
./SharedLib 
cd ../SquareEtch
echo "=== SquareEtch ==="
./SquareEtch 
cd ../VoidEtching
echo "=== VoidEtching ==="
./VoidEtching 
cd ../VolumeToLevelSets
echo "=== VolumeToLevelSets ==="
./VolumeToLevelSets 
