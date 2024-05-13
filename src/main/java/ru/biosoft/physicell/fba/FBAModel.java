package ru.biosoft.physicell.fba;

public interface FBAModel
{

    public double getFlux(String reaction);

    public void setReactionLowerBound(String rId, double lowerBound);

    public boolean getSolutionStatus();

    public double getObjectiveValue();
    
    public void runFBA();
    
    public FBAModel clone();
}