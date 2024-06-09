# Contributing to PhiMATS

## `Pointers` and `References`

- Always use `dynamic` arrays only when needed.
- Limit the usage to dynamic allocation to the interfaces with `PETSc`, i.e. within `Models` classes.
- Allocate and deallocate using `PetscMalloc1` and `PetscFree`.

## `C++` style

- Use `camelCase` for variables and `PascalCase` for function names.
- For `get_` and `set_` for getters and setters for private variables.
- Use `setAction_` for setting an action on private variable.
- For `pointers` or `references` use the style `type *` and `type &`.
