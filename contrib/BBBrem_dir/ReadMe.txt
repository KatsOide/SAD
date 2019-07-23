ReadMe.txt

Installation, should work on OS X / Unix / Linux systems

Expand  BBBrem.tgz  ( BBBrem_dir.tgz )
and keep source, for example as ~/BBBrem_dir

For gmake, make, choose a working directory for temporary output
like WORKDIR=/tmp/$LOGNAME


Example 1. Call BBBrem from another program
-------------------------------------------
for example ROOT ( https://en.wikipedia.org/wiki/ROOT )
In ROOT after loading BBRREM ( .I ~/BBBrem_dir, .L ~/BBBrem_dir/BBBrem.C ) do
double roots, k0;
BBBrem bbbrem(roots=91.2, k0=0.01);
for(auto k=0; k<10000; ++k) bbbrem.generate();
bbbrem.finish();

Example 2. BBBrem as standalone line command tool
-------------------------------------------------
cd $WORKDIR
cmake BBBrem_dir ; make
./BBBrem -h

Example 3. Use command line BBBrem to generate with binary output and read back in ROOT
---------------------------------------------------------------------------------------
BBBrem -vb 45.6 0.01 1.e5 0 6.e-11
root
.I ~/BBBrem_dir
.L ~/BBBrem_dir/BBBrem.C
BBBrem bbbrem; // default constructor
bbbrem.read_evts_binary(); // read/print first events
.q
