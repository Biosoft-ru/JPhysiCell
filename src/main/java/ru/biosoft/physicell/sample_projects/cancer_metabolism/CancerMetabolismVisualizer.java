package ru.biosoft.physicell.sample_projects.cancer_metabolism;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class CancerMetabolismVisualizer extends AgentVisualizer
{
    @Override
    public Color[] findColors(Cell cell)
    {
        if( !cell.phenotype.death.dead )
            return new Color[] {Color.black, Color.black};

        return new Color[] {Color.blue, Color.black};
    }
}