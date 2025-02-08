"""Solves the bioheat equation."""

import numpy as _np


class BioHeat:
    def __init__(self, field, limits, material, frequency, ambient_temperature=37):
        """
        Base class for solving the bioheat equation.

        Parameters
        ----------
        field: numpy.ndarray
            Acoustic pressure field in Pa.
        limits: tuple
            Contains the x,y and z limits of the simulation in the form
            (x_min, x_max, y_min, y_max, z_min, z_max).
        material: optimus.material.common.Material
            Propagation material.
        frequency: float
            Wave frequency in Hz.
        ambient_temperature: float
           Ambient temperature of the simulation in Celsius degrees.
        """

        if field.ndim != 3:
            raise TypeError("field must be 3 dimensional.")

        self.field = field
        self.material = material
        self.frequency = frequency
        self.ambient_temperature = ambient_temperature
        self.temperature = _np.zeros_like(field) + ambient_temperature
        self.amplitude = None
        self.limits = limits
        self.max_temperatures = [ambient_temperature]
        self.time_domain = _np.array([0])

        self._calculate_heat_deposition()

    def _calculate_heat_deposition(self):
        """
        Computes the the heat deposition field.

        The result is saved in self.heat_deposition in W/m^3.
        """

        self.intensity = (self.field) ** 2 / (
            2 * self.material.speed_of_sound * self.material.density
        )
        self.heat_deposition = (
            2
            * self.intensity
            * self.frequency
            * self.material.attenuation_coeff_a
            / 1000000
        )

    def _apply_finite_differences(self, delta, pulse_amplitude):
        """
        Runs one iteration of finite differences to solve the bioheat equation.
        Returns the temperature grid at the end of finite differences.

        Parameters
        ----------
        delta: tuple
            Contains the steps in the x, y and z direction and the step taken on time.
            The format is (delta_x, delta_y, delta_z, delta_time)
        pulse_amplitude: float
            Proportion of the total heat deposition given by the source.
        """

        dx, dy, dz, dt = delta
        Nz, Ny, Nx = _np.shape(self.heat_deposition)

        returning_temp = (
            _np.zeros_like(self.heat_deposition) + self.ambient_temperature
        )

        heat_constant = dt / (self.material.density * self.material.heat_capacity)
        x_constant = (
            heat_constant * self.material.thermal_conductivity * (dx ** (-2))
        )
        y_constant = (
            heat_constant * self.material.thermal_conductivity * (dy ** (-2))
        )
        z_constant = (
            heat_constant * self.material.thermal_conductivity * (dz ** (-2))
        )
        time_constant = 1 - 2 * (x_constant + y_constant + z_constant)

        returning_temp[1 : Nz - 1, 1 : Ny - 1, 1 : Nx - 1] = (
            self.heat_deposition[1 : Nz - 1, 1 : Ny - 1, 1 : Nx - 1]
            * (pulse_amplitude**2)
            * heat_constant
            + self.temperature[1 : Nz - 1, 1 : Ny - 1, 1 : Nx - 1] * time_constant
            + (
                self.temperature[1 : Nz - 1, 1 : Ny - 1, 0 : Nx - 2]
                + self.temperature[1 : Nz - 1, 1 : Ny - 1, 2:Nx]
            )
            * x_constant
            + (
                self.temperature[1 : Nz - 1, 0 : Ny - 2, 1 : Nx - 1]
                + self.temperature[1 : Nz - 1, 2:Ny, 1 : Nx - 1]
            )
            * y_constant
            + (
                self.temperature[0 : Nz - 2, 1 : Ny - 1, 1 : Nx - 1]
                + self.temperature[2:Nz, 1 : Ny - 1, 1 : Nx - 1]
            )
            * z_constant
        )

        return returning_temp

    def _select_pulse_function(
        self,
        time,
        period_duration=2e-3,
        period_repetition_interval=4e-3,
        ramp_duration=0,
        ramp_function="continuous",
    ):
        """
        Creates the pulse amplitude in the time domain in an array.

        Parameters
        ----------
        time: float
            Duration of the simulation in seconds
        period_duration: float
            Time in seconds in which the source is on.
        period_repetition_interval: float
            Duration of one pulse scheme in seconds
        ramp_duration: float
            Only used in the "turk" ramp function.
            Duration in seconds of one of the ramps portion of the pulse.
        ramp_function: str
            Function that describes the pulses.
            The options are "continuous", "turk" or "rectangular".
        """

        amplitude = _np.zeros_like(time)
        pulse_remainder = _np.remainder(time, period_repetition_interval)

        if ramp_function == "rectangular":

            indicator = pulse_remainder < period_duration
            amplitude[indicator] = 1

            return amplitude

        if ramp_function == "turk":

            indicator_1 = pulse_remainder < ramp_duration
            indicator_2 = _np.logical_and(
                pulse_remainder < period_duration - ramp_duration,
                pulse_remainder >= ramp_duration,
            )
            indicator_3 = _np.logical_and(
                pulse_remainder >= period_duration - ramp_duration,
                pulse_remainder < period_duration,
            )

            amplitude[indicator_2] = 1
            amplitude[indicator_1] = (
                -_np.cos(_np.pi * (pulse_remainder[indicator_1] / ramp_duration)) / 2
                
                + 0.5
            )
            amplitude[indicator_3] = (
                -_np.cos(
                    _np.pi
                    * (
                        (pulse_remainder[indicator_3] - period_duration)
                        / ramp_duration
                    )
                )
                / 2
                + 0.5
            )

            return amplitude

        return amplitude + 1

    def calculate_temperature(
        self,
        time,
        Nt=int(1e3),
        period_duration=2e-3,
        period_repetition_interval=4e-3,
        ramp_duration=0,
        ramp_function="continuous",
    ):
        """
        Calculates the temperature field by solving the
        bioheat equation using finite differences.
        Stores the final temperature in self.temperature.
        Stores the maximum temperature in the entire domain
        for each time iteration in self.max_temperature.

        Parameters
        ----------
        time: float
            Duration of the simulation in seconds.
        Nt: int
            Number of time steps.
        period_duration: float
            Time in seconds in which the source is on.
        period_repetition_interval: float
            Duration of one pulse scheme in seconds
        ramp_duration: float
            Only used in the "turk" ramp function.
            Duration in seconds of one of the ramps portion of the pulse.
        ramp_function: str
            Function that describes the pulses.
            The options are "continuous", "turk" or "rectangular".
        """

        xmin, xmax, ymin, ymax, zmin, zmax = self.limits
        Nz, Ny, Nx = _np.shape(self.heat_deposition)

        dx = (xmax - xmin) / Nx
        dy = (ymax - ymin) / Ny
        dz = (zmax - zmin) / Nz
        dt = time / Nt
        delta = (dx, dy, dz, dt)

        self.time_domain = _np.concatenate(
            (self.time_domain, self.time_domain[-1] + _np.array(range(Nt)) * dt)
        )

        pulse_amplitude = self._select_pulse_function(
            time=_np.array(range(Nt)) * dt,
            period_duration=period_duration,
            period_repetition_interval=period_repetition_interval,
            ramp_duration=ramp_duration,
            ramp_function=ramp_function,
        )

        for period in range(Nt):
            self.max_temperatures.append(_np.max(self.temperature))
            self.temperature = self._apply_finite_differences(
                delta, pulse_amplitude[period]
            )

    def plot_temperature(self, layer, axis, colormap_limits=None):
        """
        Plots a slice of the last temperature grid computed.

        Parameters
        ----------
        layer: int
            Layer of the slice.
        axis: int
            Axis of the layer. The options are 0, 1 or 2.
        """

        from optimus.postprocess.plot import surface_plot

        xmin, xmax, ymin, ymax, zmin, zmax = self.limits

        if axis == 0:
            plot_layer = self.temperature[layer, :, :]
            axes_labels = ("y [m]", "z [m]")
            axes_lims = (ymin, ymax, zmin, zmax)
        elif axis == 1:
            plot_layer = self.temperature[:, layer, :]
            axes_labels = ("x [m]", "z [m]")
            axes_lims = (xmin, xmax, zmin, zmax)
        elif axis == 2:
            plot_layer = self.temperature[:, :, layer]
            axes_labels = ("x [m]", "y [m]")
            axes_lims = (xmin, xmax, ymin, ymax)
        else:
            raise ValueError("Axis should be 0, 1 or 2, not " + str(axis))

        if colormap_limits == None:
            colormap_limits = [_np.min(plot_layer), _np.max(plot_layer)]

        fig = surface_plot(
            plot_layer,
            axes_lims=axes_lims,
            axes_labels=axes_labels,
            colormap="viridis",
            colormap_lims=colormap_limits,
            colorbar_label="Temperatue $[^\circ C]$",
        )

        
        
