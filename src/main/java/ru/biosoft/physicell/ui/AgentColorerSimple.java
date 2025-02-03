package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.core.Cell;

public class AgentColorerSimple implements AgentColorer
{
    private Map<Integer, Color> typeColor = new HashMap<>();
    private Map<String, Color> phaseColor = new HashMap<>();

    public void addTypeColor(Integer type, Color color)
    {
        typeColor.put( type, color );
    }

    public void addPhaseColor(String type, Color color)
    {
        phaseColor.put( type, color );
    }

    @Override
    public Color[] findColors(Cell cell)
    {
        return new Color[] {Color.black, findColor(cell)};
    }
    
    
    private Color findColor(Cell cell)
    {
        Color c = Color.white;
        int type = cell.type;
        if( typeColor.containsKey( type ) )
            c = typeColor.get( type );
        else
        {
            String phase = cell.phenotype.cycle.currentPhase().name;
            if( phaseColor.containsKey( phase ) )
                c = phaseColor.get( phase );
        }
        return c;
    }
}