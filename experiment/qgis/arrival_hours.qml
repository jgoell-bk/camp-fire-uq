<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1"
                    classificationMin="0" classificationMax="24" nodataColor="">
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0"
                         minimumValue="0" maximumValue="24">
          <item alpha="255" value="0" label="0 h" color="#cde2fb"/>
          <item alpha="255" value="4" label="4 h" color="#9ec5f4"/>
          <item alpha="255" value="8" label="8 h" color="#6da7ec"/>
          <item alpha="255" value="12" label="12 h" color="#3987e5"/>
          <item alpha="255" value="16" label="16 h" color="#256abf"/>
          <item alpha="255" value="20" label="20 h" color="#184f95"/>
          <item alpha="255" value="24" label="24 h" color="#0d366b"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
  </pipe>
</qgis>
