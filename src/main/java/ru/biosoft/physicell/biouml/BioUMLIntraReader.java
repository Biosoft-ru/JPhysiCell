package ru.biosoft.physicell.biouml;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Map;

import org.w3c.dom.Element;

import biouml.model.Diagram;
import biouml.plugins.sbml.SbmlModelFactory;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.xml.IntracellularReader;
import ru.biosoft.physicell.xml.ModelReaderSupport;

public class BioUMLIntraReader extends ModelReaderSupport implements IntracellularReader
{
    @Override
    public void readIntracellular(Path path, Element el, Model model, CellDefinition cd) throws Exception
    {
        Phenotype p = cd.phenotype;
        String type = getAttr( el, "type" );
        if( type.equals( "dfba" ) )
        {
            return;
        }
        IntracellularODEBioUML intracellular = new IntracellularODEBioUML( model, cd );
        p.intracellular = intracellular;
        for( Element child : getAllElements( el ) )
        {
            String tag = child.getTagName();
            if( tag.equals( "intracellular_dt" ) )
            {
                double dt = getDoubleVal( child );
                intracellular.setDT( dt );
            }
            else if( tag.equals( "map" ) )
            {
                String species = getAttr( child, "sbml_species" );

                if( hasAttr( child, "PC_substrate" ) )
                {
                    String substrate = getAttr( child, "PC_substrate" );
                    intracellular.addPhenotypeSpecies( substrate, species );
                }
                else if( hasAttr( child, "PC_phenotype" ) )
                {
                    String code = getAttr( child, "PC_phenotype" );
                    intracellular.addPhenotypeSpecies( code, species );
                }
            }
            else if( tag.equals( "sbml_filename" ) )
            {
                String name = getVal( child );
//                if (name.contains("/"))
//                    name = name.substring( name.lastIndexOf( "/" ) );
                InputStream stream = model.getClass().getResourceAsStream( name );
                intracellular.setDiagram( readSBML( stream, name ) );
            }
        }
    }

    public Diagram readSBML(InputStream stream, String name) throws Exception
    {
        return SbmlModelFactory.readDiagram( stream, name, null, name );
    }

    @Override
    public void setAdditionalFiles(Map<String, File> additional)
    {
    }
}