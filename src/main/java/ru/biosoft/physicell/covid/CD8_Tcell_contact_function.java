package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.Contact;
import ru.biosoft.physicell.core.standard.StandardModels;

public class CD8_Tcell_contact_function extends Contact
{

    @Override
    public void execute(Cell pC1, Phenotype p1, Cell pC2, Phenotype p2, double dt)
    {
        // std::cout << pC1 << " " << pC1.type_name 
        // << " contact with " << pC2 << " " << pC2.type_name << std::endl; 
        // elastic adhesions 
        StandardModels.standard_elastic_contact_function( pC1, p1, pC2, p2, dt );

        // increase contact time of cell you are attacking 
        //        #pragma omp critical
        //            {
        pC2.customData.set( "TCell_contact_time", pC2.customData.get( "TCell_contact_time" ) + dt );
        //            }
    }
}