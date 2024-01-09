package ru.biosoft.physicell.core;

@FunctionalInterface
public interface PhaseExit
{
    public void execute(Cell pCell, Phenotype phenotype, double dt);
}
