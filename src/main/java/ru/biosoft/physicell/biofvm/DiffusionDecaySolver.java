package ru.biosoft.physicell.biofvm;

public abstract class DiffusionDecaySolver
{
    public abstract String getName();
    public abstract void solve(Microenvironment m, double dt) throws Exception;
}
