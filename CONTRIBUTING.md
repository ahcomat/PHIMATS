# Contributing to PhiMATS

## `C++` style

- Always use `dynamic` arrays only when needed.
- Limit the usage to dynamic allocation to the interfaces with `PETSc`, i.e. within `Models` classes.
- Allocate and deallocate using `PetscMalloc1` and `PetscFree`.
- For `pointers` or `references` use the style `type *` and `type &`.
