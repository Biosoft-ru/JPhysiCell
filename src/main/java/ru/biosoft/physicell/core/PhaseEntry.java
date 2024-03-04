package ru.biosoft.physicell.core;

@FunctionalInterface
public interface PhaseEntry
{
    public void execute(Cell pCell, Phenotype phenotype, double dt);
}