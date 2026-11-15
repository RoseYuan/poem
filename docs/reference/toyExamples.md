# Toy embedding examples

Toy example 2D embeddings of elements of different classes, with varying
mixing and spread. Graphs 1-3 all have 20 elements of each of 4 classes,
but that are mixed in different fashion in the embedding space. Graphs
4-7 all have 100 elements of class1 and 60 of class2, and the class1
elements vary in their spread. The code to generate the dataset is at
system.file('inst/scripts/', 'graph_example.R', package='poem').

## Usage

``` r
data(toyExamples)
```

## Format

### `toyExamples`

A data frame.

- graph:

  The name of the embedding to which the element belongs.

- x, y:

  Coordinates in the 2D embedding.

- class:

  The class to which the element belongs.
