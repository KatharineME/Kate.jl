Julia and I met in 2018. We started dating recently. More I spend time with her, more things I discover.

## Not-so-obvious concepts from https://juliadocs.github.io/Julia-Cheat-Sheet

Swap pointers

``` julia
x, y = y, x
```

Chain from right to left

``` julia
x = y = z = 1

z = 1; y = z; x = y
```

Divide in reverse

``` julia
7\3 == 3/7
```

Think of REPL to have 3 extra modes: `;` , `]` , and `?` .

Preallocate a Float64 array of size N

``` julia
array = Float64[]

sizehint!(array, N)
```

Make a Float64 array of size N with undef

``` julia
Array{Float64, 1}(undef, N)
```

All function arguments are passed by pointers; messing with an argument value within a function messes the value up everywhere else.

 `struct` = `Type` = `class` - functions - mutable (unless `mutable struct` ) .

Use immutable types for threaded applications to avoid synchronizing the objects.

Types can be parametrized; `T1` holding child type of `T2` is `T1{T<:T2}` .

Only `Tuple` is covariate: `Tuple{Float64} <: Tuple{Real}` returns `true` .

There is only 1 `missing` and 1 `nothing` ( `missing === missing` and `nothing === nothing` both return `true` ).

 `missing` is contageous; any evaluation of it returns `missing` .

Filter missings

``` julia
 collect(skipmissing(collection))
```

Get exported names

``` julia
names(ModuleName)
```

Make `macro` like making `function` .

Use `macro` like using `function` or `@macroname expression` .

Built-in macros

```julia
@test
@test_approx_eq
@assert
@which
@time
@elapsed
@allocated
@profile
@spawn
@spawnat
@async
@distibuted
@everywhere
```julia

## map vs collection comprehension

 `map` preserves the type of a collection while list comprehension does not.

``` julia
set = Set[1, 2, 3]

map(x -> x + 1, set) # returns a Set

[x + 1 for x in set] # returns an Array
```

 `map` is dispatch friendly.

## ===, isequal, ==

Everything lives in a physical location in a computer. Here is a string "Kwat". Here is another string "Kwat". Although these two strings have the same values (text "Kwat"), they are located in different locations in a computer.
 `===` checks if two things have the same physical locations, while `==` checks if two things have the same values.

``` julia
 pointer1 = "Kwat"

 pointer2 = "Kwat"

 pointer1 === pointer2 # false

 pointer1 == pointer2 # true
 ```

 `isequal` is very much like `==` , but with some differences.

``` julia
0.0 == -0.0 # true

isequal(0.0, -0.0) # false

NaN == NaN # false

isequal(NaN, NaN) # true
```

These differences make `isequal` more suited for sorting because `-0.0` s, `0.0` s, and `NaN` s will stick together.

## Performance tips

Avoid global variables.

Avoid abstract types in collections.

Type explicitly and specifically.

Use immutable types as much.

Preallocate large arrays.

Preallocate data structure for results.

Free up memory by setting a pointer to `nothing` .

Access array along columns.

Mutate as much!

Broadcast instead of using collection comprehensions.

Avoid `try-catch` .

Aovid string interpolation in I/O.

Vectorizing is often slower (if I'm good).

