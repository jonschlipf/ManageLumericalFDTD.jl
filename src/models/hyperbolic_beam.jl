export prepare_hyperbolic_beam_2D
using PyCall
function prepare_hyperbolic_beam_2D(params,filename="test.fsp")
    sys=pyimport("sys")
    os=pyimport("os")
    sys.path=cat(sys.path,ENV["LUMERICAL_PYTHON_PATH"],dims=1)
    lu=pyimport("lumapi")
    fdtd=lu.FDTD(hide=true)#hide=params[:hide])
    zspace=params["λmax"]*5
    xspace=params["waist"]
    angle_in=asind(3.5*asin(params["θ"]))
    zsize=zspace*3+params["tsub"]+params["thmm"]+params["tvs"]
    x_min=-xspace*1e-9,
    x_max=xspace*1e-9+params["tsub"]*1.2*atand(angle_in)*1e-9,
     y_max=2*zspace*1e-9,
     y_min=-1e-9*(zspace+params["tsub"]+params["thmm"]+params["tvs"]),
    fdtd.addfdtd(dimension=("2D"),
                 x_min=x_min,
                 x_max=x_max,
                 y_min=y_min,
                 y_max=y_max,
                 dt_stability_factor=params["dt"],
                 auto_shutoff_min=params["aso"],
                 simulation_time=params["t"],
                 x_min_bc="PML", 
                 y_min_bc="PML", 
                 x_max_bc="PML", 
                 y_max_bc="PML", 
                 mesh_accuracy=params["ma"],
                 )
    fdtd.addgaussian(injection_axis="y",
                     direction="Backward",
                     x=0,
                     y=zspace*1e-9,
                     x_span=2*(xspace*1e-9+params["tsub"]*1.2*atand(angle_in)*1e-9),
                     angle_theta=-1*params["θ"],
                     polarization_angle=params["pol"],
                     waist_radius_w0=params["waist"]*1e-9*.5,
                     frequency_dependent_profile=true,
                     override_global_source_settings=false,
                     )
    fdtd.setglobalsource("wavelength start",params[:λmin]*1e-9)
    fdtd.setglobalsource("wavelength stop",params[:λmax]*1e-9)
    fdtd.setglobalmonitor("frequency points",params[:fpoints])
    fdtd.setglobalmonitor("use source limits",true)
    fdtd.setglobalmonitor("use wavelength spacing",true)
    fdtd.addpower(dimension=("2D"),
                  name="field",
                 x_min=x_min,
                 x_max=x_max,
                 y_min=y_min,
                 y_max=y_max,
                  )
    fdtd.addrect(material="Si (Silicon) - Palik",
                 x_min=x_min,
                 x_max=x_max,
                 y_max=-(params["tvs"]+params["thmm"])*1e-9,
                 y_max=-(params["tvs"]+params["thmm"]+params["sub"])*1e-9
                 )
    fdtd.addrect(material="Ge (Germanium) - CRC",
                 x_min=x_min,
                 x_max=x_max,
                 y_max=-(params["thmm"])*1e-9,
                 y_max=-(params["tvs"]+params["thmm"])*1e-9
                 )
    fdtd.addrect(material="Ge (Germanium) - CRC",
                 x_min=x_min,
                 x_max=x_max,
                 y_max=0,
                 y_max=-params["thmm"]*1e-9
                 )
    fdtd.save(filename)
    fdtd.close()
end
