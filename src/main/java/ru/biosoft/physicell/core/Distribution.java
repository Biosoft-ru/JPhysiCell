package ru.biosoft.physicell.core;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

public class Distribution
{

    private String name;
    private Map<String, Double> parameters = new HashMap<>();

    public Distribution(String name)
    {
        this.name = name;
    }
    
    public String getName()
    {
        return name;
    }

    public void addParameter(String parameter, Double value)
    {
        parameters.put( parameter, value );
    }
    
    public double getValue(String parameter)
    {
        return parameters.get( parameter );
    }
    
    public Distribution clone()
    {
        Distribution result = new Distribution(name);
        for (Entry<String, Double> e: parameters.entrySet())
            result.addParameter( e.getKey(), e.getValue() );
        return result;
    }
}