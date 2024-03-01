package ru.biosoft.physicell.ode;

import java.util.HashMap;
import java.util.Map;

import biouml.plugins.simulation.java.JavaBaseModel;

public class ODEModel extends JavaBaseModel
{
    protected Map<String, Integer> variableIndex = new HashMap<>();

    protected void initMap()
    {

    }

    public double getCurrentValue(String name) throws Exception
    {
        int index = variableIndex.get( name );
        return getCurrentValues()[index];
    }

    public void setCurrentValue(String name, double value) throws Exception
    {
        int index = variableIndex.get( name );
        double[] vals = getCurrentValues();
        vals[index] = value;
        setCurrentValues( vals );
    }

    public ODEModel clone()
    {
        ODEModel result = (ODEModel)super.clone();
        result.variableIndex = new HashMap<String, Integer>();
        result.variableIndex.putAll( variableIndex );
        return result;
    }
}
