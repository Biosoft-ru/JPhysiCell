package ru.biosoft.biofvm.cell;

@FunctionalInterface
public interface PhaseExit
{
    public void execute(Cell pCell, Phenotype phenotype, double dt);
}
