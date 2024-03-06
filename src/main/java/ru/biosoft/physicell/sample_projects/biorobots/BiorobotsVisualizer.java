package ru.biosoft.physicell.sample_projects.biorobots;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.AgentVisualizer;
import ru.biosoft.physicell.xml.ModelReader;

public class BiorobotsVisualizer extends AgentVisualizer
{
    Map<Integer, Color> colors = new HashMap<>();

    public BiorobotsVisualizer(Model model)
    {
        int cargoType = CellDefinition.getCellDefinition( "cargo cell" ).type;
        int workerType = CellDefinition.getCellDefinition( "worker cell" ).type;
        int directorType = CellDefinition.getCellDefinition( "director cell" ).type;
        colors.put( cargoType, ModelReader.readColor( model.getParameter( "cargo_color" ) ) );
        colors.put( workerType, ModelReader.readColor( model.getParameter( "worker_color" ) ) );
        colors.put( directorType, ModelReader.readColor( model.getParameter( "director_color" ) ) );
    }

    @Override
    public Color[] findColors(Cell cell)
    {
        return new Color[] {colors.get( cell.type )};
    }
}