package ru.biosoft.physicell.core;

import java.util.HashMap;
import java.util.Map;

public class CellContainerRegistry
{
    private static CellContainer[] containers = new CellContainer[] {new CellContainer(), new CellContainerExperimental(),
            new CellContainerParallel()};

    private static Map<String, Class<? extends CellContainer>> containersMap = new HashMap<>();

    private static boolean isInit = false;

    private static void init()
    {
        containersMap = new HashMap();
        for( CellContainer container : containers )
            containersMap.put( container.getName(), container.getClass() );
        isInit = true;
    }

    public static String[] getAvailableContainers()
    {
        return containersMap.keySet().toArray( new String[containersMap.size()] );
    }

    public static Class<? extends CellContainer> getClass(String name) throws IllegalArgumentException
    {
        if( !isInit )
            init();
        
        if( !containersMap.containsKey( name ) )
            throw new IllegalArgumentException( "Can not find cell container " + name );
        return containersMap.get( name );
    }

    public static CellContainer createCellContainer(String name) throws Exception
    {
        Class<? extends CellContainer> clazz = getClass( name );
        return clazz.getConstructor( null ).newInstance( null );
    }
}
