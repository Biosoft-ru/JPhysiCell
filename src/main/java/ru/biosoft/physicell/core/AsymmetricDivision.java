package ru.biosoft.physicell.core;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.biofvm.VectorUtil;

public class AsymmetricDivision implements Cloneable
{
    private Map<String, Integer> cellTypes = new HashMap<>();
    private boolean enabled = false;
    private double[] probabilities = new double[] {0.0};
    private String[] names = new String[] {"defaults"};
    
    public void initialize(Model model)
    {
        int size = model.getDefinitionsCount();
        for( int i = 0; i < size; i++ )
            cellTypes.put(  model.getCellDefinition( i ).name, i );
        probabilities = VectorUtil.resize( probabilities, size );
        names = new String[size];
        for (int i=0; i<size; i ++)
            names[i] = model.getCellDefinition( i ).name;
    }
    
    public double getTotalProbability()
    {
        return Arrays.stream( probabilities ).sum();
    }

    public double getProbability(String name)
    {
        Integer index = cellTypes.get( name );
        if( index == null )
            throw new IllegalArgumentException( "Incorrect cell type " + name );
        return probabilities[index];
    }
    
    public double getProbability(int index)
    {
        return probabilities[index];
    }

    public void setProbability(int index, double val)
    {
        probabilities[index] = val;
    }
    
    public void setProbability(String name, double val)
    {
        Integer index = cellTypes.get( name );
        if( index == null )
            throw new IllegalArgumentException( "Incorrect cell type " + name );

        probabilities[index] = val;
    }

    public double[] getProbabilities()
    {
        return probabilities;
    }
    
    public String[] getNames()
    {
        return names;
    }
    
    @Override
    public AsymmetricDivision clone()
    {
        try
        {
            AsymmetricDivision result = (AsymmetricDivision)super.clone();
            result.probabilities = this.probabilities.clone();
            result.names = names;
            result.cellTypes = new HashMap<>(cellTypes);
            result.enabled = enabled;
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
    
    public void setEnabled(boolean enabled)
    {
        this.enabled = enabled;
    }
    
    public boolean isEnabled()
    {
        return enabled;
    }
}