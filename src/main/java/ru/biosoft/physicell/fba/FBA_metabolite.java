package ru.biosoft.physicell.fba;

public class FBA_metabolite
{
    String id;
    String name;

    FBA_metabolite(String id)
    {
        this.id = id;
    }

    FBA_metabolite()
    {

    }

    String getId()
    {
        return this.id;
    }

    void setName(String value)
    {
        this.name = value;
    }

    String getName()
    {
        return this.name;
    }
}
