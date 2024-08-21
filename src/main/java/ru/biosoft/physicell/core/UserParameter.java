package ru.biosoft.physicell.core;

public class UserParameter
{
    private String name;
    private String description;
    private String value;
    
    public UserParameter(String name, String description, String value)
    {
        this.name  = name;
        this.description = description;
        this.value = value;
    }
    
    public String getName()
    {
        return name;
    }
    
    public String getValue()
    {
        return value;
    }
   
    public String getDescription()
    {
        return description;
    }
    
}
