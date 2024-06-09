# GOAL

## Directory Structure

- **GOAL/**
  C++ source code for the core algorihtm implementations and interfaces.

- **Test/**
  entry of console programs which invoke the core algorithm API for batch black-box testing.
  it can also be regarded as the sample code.

- **Data/**
  deploy directory for executables.

  - **OCM/**
    environment setting, configurations, and data for testing each problem.

    - **Solution/**
    save running logs and output files of this problem.

    - **Visualization/**
    save visualization files of this problem.
    
    - **./0.InstanceList.txt**
      the list of instances to test.

    - **./0.Baseline.txt**
      the best known results for each instance.

    - **./SomeInstance.txt**
      instance file.

- **Doc/**
  documents for development.

- **Build/**
  project configurations for build toolsets.


## Build Instruction

### Linux

The solver is developped with g++ v11.4.0,  cmake v2.8 in the Linux system.

To use it, you only need the C++ standard library and a g++ compiler.

To build, use:

```
cd Build/Linux
mkdir build
cd build
cmake ../
make
```

To run, use:

```
./Build/Linux/build/OCM
```

The problem is piped in via STDIN.

The solution is emitted via STDOUT.

## Contact

Kong Qi(k969774646@gmail.com)

