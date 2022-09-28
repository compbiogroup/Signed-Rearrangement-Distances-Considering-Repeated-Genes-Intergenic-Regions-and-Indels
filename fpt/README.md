# Fixed Parameter algorithm for Signed Minimum Common Intergenic String Partition

Two fixed parameter algorithms for the Signed Minimum Common Intergenic String Partition problem.

## Usage

Compile the code by running `make` and run with: 
```bash
      java --class-path <bin folder path> main.MCSP <input_path> <output_path> <NC|CC> 
```

The NC option will run the partition cost used in FPT1 and CC will use the partition cost used in FPT2, according to the description in the paper.
