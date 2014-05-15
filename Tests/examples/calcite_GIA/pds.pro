<?xml version="1.0"?>
<ogs6>
<processParameters>
	<GROUNDWATER_FLOW TimeGroupID="0"/>
	<HEAD_TO_ELEMENT_VELOCITY TimeGroupID="0"/>
    <REACT_GIA TimeGroupID="0"/>
</processParameters>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="15" epsilon="1e-4">
        <out>HEAD</out>
        <out>VELOCITY</out>
        <out>Tracer</out>
        <problems>
            <M name="GROUNDWATER_FLOW" type="GROUNDWATER_FLOW">
                <out>HEAD</out>
            </M>
            <M name="HEAD_TO_ELEMENT_VELOCITY" type="HEAD_TO_ELEMENT_VELOCITY">
                <in>HEAD</in>
                <out>VELOCITY</out>
            </M>
            <M name="REACT_GIA" type="REACT_GIA">
                <in>VELOCITY</in>
                <out>Tracer</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

