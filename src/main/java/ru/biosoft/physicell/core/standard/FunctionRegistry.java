package ru.biosoft.physicell.core.standard;

import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.core.CellFunctions.CellDivision;
import ru.biosoft.physicell.core.CellFunctions.Contact;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.CellFunctions.DistanceCalculator;
import ru.biosoft.physicell.core.CellFunctions.Function;
import ru.biosoft.physicell.core.CellFunctions.MembraneInteractions;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdateVelocity;
import ru.biosoft.physicell.core.CellFunctions.VolumeUpdate;
import ru.biosoft.physicell.core.CellFunctions.Instantiator;
import ru.biosoft.physicell.core.CellFunctions.set_orientation;
import ru.biosoft.physicell.sample_projects.pred_prey_farmer.AvoidBoundariesRule;
import ru.biosoft.physicell.sample_projects.pred_prey_farmer.WrapBoundariesRule;

public class FunctionRegistry
{
    private static Map<String, Function> mapping;
    private static boolean isInit = false;

    public static Function getFunction(String name)
    {
        if( !isInit )
            initMapping();
        return mapping.get( name );
    }

    public static void initMapping()
    {
        mapping = new HashMap<>();
        Function[] allFunctions = getAllFunctions();
        for( Function f : allFunctions )
            mapping.put( f.getName(), f );
        isInit = true;
    }

    public static Function[] getAllFunctions()
    {
        return new Function[] {new StandardVolumeUpdate(), new O2based(), new AvoidBoundariesRule(), new WrapBoundariesRule(),
                new StandardUpdateVelocity(), new StandardElasticContact(), new StandardElasticContact(), new UpOrientation(),
                new DomainEdgeAvoidance(), new DomainEdgeDistance(), new Chemotaxis(), new AdvancedChemotaxis(),
                new AdvancedChemotaxisNormalized(), new StandardAssymetricDivision()};
    }

    public static VolumeUpdate[] getVolumeFunctions()
    {
        return new VolumeUpdate[] {new StandardVolumeUpdate()};
    }

    public static UpdatePhenotype[] getUpdatePhenotypeFunctions()
    {
        return new UpdatePhenotype[] {new O2based()};
    }

    public static CustomCellRule[] getCustomRules()
    {
        return new CustomCellRule[] {new AvoidBoundariesRule(), new WrapBoundariesRule()};
    }

    public static UpdateMigrationBias[] getUpdateMigrationFunctions()
    {
        return new UpdateMigrationBias[] {new Chemotaxis(), new AdvancedChemotaxis(), new AdvancedChemotaxisNormalized()};
    }

    public static UpdateVelocity[] getVelocityFunctions()
    {
        return new UpdateVelocity[] {new StandardUpdateVelocity()};
    }

    public static Contact[] getContactFunctions()
    {
        return new Contact[] {new StandardElasticContact()};
    }

    public static Instantiator[] getIntsnatiateFunctions()
    {
        return new Instantiator[] {};
    }

    public static set_orientation[] getOrientationFunctions()
    {
        return new set_orientation[] {new UpOrientation()};
    }

    public static MembraneInteractions[] getMembraneInteractionFunctions()
    {
        return new MembraneInteractions[] {new DomainEdgeAvoidance()};
    }

    public static DistanceCalculator[] getDistanceCalculatorFunctions()
    {
        return new DistanceCalculator[] {new DomainEdgeDistance()};
    }
    
    public static CellDivision[] getDivisionFunctions()
    {
        return new CellDivision[] {new StandardAssymetricDivision()};
    }
}
