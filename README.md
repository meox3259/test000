Detect marcrosatellite test code

```bash
rm -rf build && mkdir build && cd build && cmake .. && make && cd ..
```

run
```bash
./build/bin/suffix -v -f data/test2.fa -r result.txt
```