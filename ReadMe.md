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

### Windows

- Install Protobuf.
  - Set environment variable `PROTOBUF_HOME` to Protobuf installation path (where there are `bin/`, `lib/`, `include/` directories).
  - Append `%PROTOBUF_HOME%` to the envrionment variable `Path`.
  - Run `/Plugin/Protobuf/0.Install.*.bat`.
- Install Zlib.
  - Set environment variable `ZLIB_HOME` to Zlib installation path (where there are `bin/`, `lib/`, `include/` directories).
  - Append `%ZLIB_HOME%` to the envrionment variable `Path`.
- Install Gurobi.
  - Set environment variable `GUROBI_HOME` to Gurobi installation path (where there are `bin/`, `lib/`, `include/` directories).
    - Already done by the Gurobi installer.
  - Append `%GUROBI_HOME%` to the envrionment variable `Path`.
    - Already done by the Gurobi installer.
- Install core API.
  - Add all files under `/GOAL/` into your project.
  - Add the path where you place `/GOAL/` to your project include directory.
  - Include the `/GOAL/Optimization/*/Problem.h` in your code.
  - Refer to `/Test/Test*.cpp` for sample usage.
