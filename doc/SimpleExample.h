/*!
\page SimpleExample A Simple Example
\dontinclude examples/simple-example.cpp
\section A Simple Example
\c libfbi was designed to be straightforward and easy to use. In brief, it
requires the user to specify a simple adaption layer that is used to call one of
two static workhorse functions.

First we have to include libfbi
\skip #include
\until intersect 

Assume the data consists of a set of structs of type
\skip struct
\until };

A possible suitable container would be
\skip std
\until locations;


Now, assume you would like to know which locations are in walking distance of
each other, and that you are willing to walk at most a mile (city block
distance). Hence, you draw a rectangle with one mile sides and center it over
each location. The task now is to determine for which locations the rectangles 
overlap, indicating that they are in walking distance.from each other.

Let's define a box operator for the locations. First, we need to define some
technical traits (to define the comparison operators, etc). As \c Location only
holds doubles, this is straightforward: just use the \c TraitsGenerator:

\skip namespace
\until end namespace



Now we can define a box generator for a \c Location object:
\skip struct
\until };

The template argument \c N defines the box dimension and we use 0 for x and 1
for y. Specifying the template functions, we get
\skip template
\until }

for the x dimension and

\skip template
\until }

for the y dimension. Determining the adjacency list that holds the indexes of
all elements whose boxes intersect is then a simple matter of calling the single
set intersection function

\skip main
\until }

The resulting adjecency list can subsequently be used to calucate connected
components or similar.
*/
