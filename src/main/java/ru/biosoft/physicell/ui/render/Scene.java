package ru.biosoft.physicell.ui.render;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Scene
{
    private Map<Integer, List<Mesh>> layers = new HashMap<>();

    private List<Mesh> spheres = new ArrayList<Mesh>();

    public void addDisk(Mesh mesh, int layerCode)
    {
        layers.computeIfAbsent( layerCode, i -> new ArrayList<>() ).add( mesh );
    }

    public void addSphere(Mesh mesh)
    {
        spheres.add( mesh );
    }

    public void clearLayer(int layerCode)
    {
        List<Mesh> layer = layers.get( layerCode );
        if( layer != null )
            layer.clear();
    }
    
    public void clear()
    {
        spheres.clear();
        for( List<Mesh> layer : layers.values() )
            layer.clear();
        layers.clear();
    }

    public int getSpheresCount()
    {
        return spheres.size();
    }

    public Iterable<Mesh> getSpheres()
    {
        return spheres;
    }

    public Iterable<Mesh> getLayer(int layerCode)
    {
        return layers.getOrDefault( layerCode , new ArrayList<>());
    }

    
    public void sortSpheres()
    {
        spheres.sort( null );
    }
}
