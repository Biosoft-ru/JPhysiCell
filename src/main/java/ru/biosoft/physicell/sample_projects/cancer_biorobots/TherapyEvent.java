package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Model.Event;

public class TherapyEvent extends Event
{
    public TherapyEvent(double executionTime)
    {
        super( executionTime );
    }

    @Override
    public void execute(Model model) throws Exception
    {
        System.out.println( "Therapy started!" );
        model.setSaveInterval( model.getParameterDouble( "save_interval_after_therapy_start" ) ); // 3.0; 
        CancerBiorobots.introduce_biorobots( model );
    }
}