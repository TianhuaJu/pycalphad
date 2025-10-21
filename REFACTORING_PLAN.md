# Refactoring Plan for UEM Implementation

## Issues Identified by Project Lead

1. **Implementation is too verbose** - Should be simpler with less code
2. **Tests not following conventions** - Need to be in-tree, use pytest fixtures
3. **Missing numeric verification** - Need analytic solutions or code-to-code comparison
4. **Documentation format wrong** - Need RestructuredText (.rst), not Markdown
5. **Documentation too verbose** - Should be more concise
6. **Examples format** - Should be Jupyter notebooks

7. **References are WRONG** - Author is **Tianhua Ju**, not Chou K.C.

## Correct References

Based on research, the correct author attribution is:

**Lead Author**: Tianhua Ju
**Co-authors**: Xueyong Ding, Long Zhang, Xinlin Yan, and others

**Publications (2020)**:
1. Ju, T., Ding, X., Zhang, L., Chen, W., Wang, B., & Yan, X. (2020). "On the definition of the components' difference in properties in the unified extrapolation model." *Fluid Phase Equilibria*.

2. Ju, T., Ding, X., Yan, X., Liu, C., Zhang, X., & Zhang, L. (2020). "New expression for property difference in components for the Unified Extrapolation Model." *Journal of Molecular Liquids*.

## Refactoring Tasks

### Priority 1: Core Code Simplification
- [ ] Simplify uem_symbolic.py - remove excessive comments
- [ ] Keep only essential docstrings
- [ ] Reduce code duplication
- [ ] Target: 50% reduction in code volume

### Priority 2: Correct References
- [ ] Update all citations to Tianhua Ju (2020)
- [ ] Remove all references to "Chou K.C."
- [ ] Add proper DOI links if available

### Priority 3: Test Restructuring
- [ ] Delete test_uem_validation.py from root
- [ ] Move all tests to pycalphad/tests/test_model.py
- [ ] Use @select_database decorator
- [ ] Add numeric verification tests
- [ ] Verify against analytic solutions for simple cases

### Priority 4: Documentation Format
- [ ] Convert Markdown docs to RestructuredText (.rst)
- [ ] Reduce verbosity by 70%
- [ ] Focus on essential information only
- [ ] Follow Sphinx conventions

### Priority 5: Examples
- [ ] Create 1-2 Jupyter notebooks
- [ ] Show UEM vs Muggianu comparison
- [ ] Include plots and visualizations
- [ ] Reference the Ju et al. papers

### Priority 6: Integration
- [ ] Consider integrating UEM into base Model class
- [ ] Similar to how Kohler-Toop is integrated
- [ ] Make it selectable via parameter

## Target Metrics

- **Code reduction**: From ~1,500 lines to ~300-500 lines
- **Doc reduction**: From ~2,500 lines to ~500 lines
- **Tests**: All in-tree, following pytest conventions
- **Examples**: 2 Jupyter notebooks instead of long .py files

## Timeline

1. Fix references (immediate)
2. Simplify code (1-2 hours)
3. Restructure tests (1 hour)
4. Convert documentation (1 hour)
5. Create Jupyter notebooks (1 hour)
6. Final review and commit (30 min)

Total estimated: 4-6 hours of focused work

## Next Steps

Start with fixing all the reference citations, then systematically work through the priorities.
