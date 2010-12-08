/*!
\page AdvancedExample A more advanced example
\dontinclude examples/soccer-example.cpp

This toy example illustrates the use of two sets in the fast box intersection
procedure, where the underlying sets have different types. The sample use case
mimicks the task of a soccer scout, searching soccer players that fit a
position, age and salary description that clubs are interested in and, vice
versa, also determines the set of clubs each player could play in, given his
player profile.

First we set up the includes
\skip #include
\until intersect 

A soccer player is represented by the following struct
\skip struct
\until };
whereas clubs are modeled using 
\skip struct
\until };

Now we need to define the respective box operators, that return intervals in
each dimension. As a precondition, we need to define the traits for each of the
classes:

\skip namespace
\until end namespace

The box generator class definition for the \c Player class is

\skip struct
\until };

and we need to specialize the template member function for each dimension.
For the player, we choose a dimension order of \c position, \c price, and \c
age. This yields

\skip template
\until }

for the position,

\skip template
\until }

for the price, and

\skip template
\until }

for the age dimension.

We repeat the procedure for the \c Club class, however, for illustration
purposes, we flip the order of the dimension (and unflip it in the call to
\c intersect later). The definition of the box operator is

\skip class
\until };

and the necessary template specifications are

\skip template
\until }

for the position,

\skip template
\until }

for the age, and

\skip template
\until }

for the price dimension (note the different dimension order!).

To illustrate our type definitions, we write a short \c main program

\skip main
\until clubs;

We now make use of one convenient property of the libfbi API, which allows us to
pass multiple box operators into the intersection procedure. We instanciate \c
Club boxes for every position for which players are sought.

\skip std
\until 11, std::make_pair(19, 24)));

Finding suitable players for all clubs and suitable clubs for all players now is
a simple matter of calling the nested \c intersect function. Note the flipped
dimension indexes in the inner \c SetB template specification.

\skip auto
\until );

That's it!

\skip return
\until }

*/
