package ru.biosoft.physicell.ui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DensityState
{
    private String name;
    private List<double[]> densities = new ArrayList<>();
    private Map<String, Integer> names = new HashMap<String, Integer>();

    public DensityState(String name)
    {
        this.name = name;
    }
    public void addDensity(String name, double[] density)
    {
        int index = densities.size();
        densities.add( density );
        names.put( name, index );
    }

    public double[] getDensity(String name)
    {
        Integer index = names.get( name );
        if( index == null )
            return null;
        return densities.get( index );
    }
    
    public String[] getSubstrates()
    {
        return names.keySet().toArray( String[]::new );
    }
    
    public String getName()
    {
        return name;
    }
}