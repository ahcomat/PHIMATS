# Contributing to PhiMATS

## General practices

- For `std::vector`, unless extremely necessary, don't use `pushback()` and instead use `resize()`. This avoids changing the size of the vector by mistake.


## `Pointers` and `References`

- Always use `dynamic` arrays only when needed.
- Limit the usage to dynamic allocation to the interfaces with `PETSc`, i.e. within `Models` classes.
- Allocate and deallocate using `PetscMalloc1` and `PetscFree`.

## `C++` style

- For includes, first add standard includes, then custom ones. Add a blank line after them.
- Element specific data begin with `nEl`.
- Total model data begin with `nTot`.
- Element set data begin with `n`.
- Use `camelCase` for variables and `PascalCase` for function names.
- For `get_` and `set_` for getters and setters for private variables.
- Use `setAction_` for setting an action on private variable.
- For `pointers` or `references` use the style `type *` and `type &`.
