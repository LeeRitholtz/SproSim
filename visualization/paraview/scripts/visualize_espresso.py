#!/usr/bin/env python3
"""
ParaView Python script for automated espresso extraction visualization.

This script loads SproSim espresso simulation data and creates a comprehensive
visualization showing particle extraction states and flow fields over time.

Usage:
    1. Open ParaView
    2. Go to Tools -> Python Shell
    3. Run: exec(open('visualize_espresso.py').read())

Or from command line:
    pvpython visualize_espresso.py /path/to/brewing_simulation.pvd
"""

import sys
import os
from paraview.simple import *

def setup_espresso_visualization(pvd_file_path, output_dir=None):
    """
    Set up comprehensive espresso extraction visualization.

    Args:
        pvd_file_path: Path to the brewing_simulation.pvd file
        output_dir: Optional directory to save screenshots/animations
    """

    # Clear any existing pipeline
    Delete(GetSources())

    # Load the time series data
    print(f"Loading data from: {pvd_file_path}")
    reader = PVDReader(FileName=pvd_file_path)

    # Get data information
    reader.UpdatePipeline()
    data_info = reader.GetDataInformation()
    point_data = data_info.GetPointDataInformation()

    print(f"Number of timesteps: {len(reader.TimestepValues)}")
    print(f"Number of points per timestep: {data_info.GetNumberOfPoints()}")
    print("Available data arrays:")
    for i in range(point_data.GetNumberOfArrays()):
        array_info = point_data.GetArray(i)
        print(f"  - {array_info.Name} ({array_info.GetNumberOfComponents()} components)")

    # Create particle visualization with spheres
    print("Creating particle glyphs...")
    sphere_source = Sphere(Radius=0.0002, ThetaResolution=8, PhiResolution=6)
    particle_glyphs = Glyph(Input=reader, GlyphType=sphere_source)
    particle_glyphs.ScaleMode = 'scalar'
    particle_glyphs.ScaleFactor = 5.0  # Make particles more visible
    particle_glyphs.GlyphMode = 'All Points'
    particle_glyphs.MaximumNumberOfSamplePoints = 5000

    # Create render view
    render_view = CreateView('RenderView')
    render_view.ViewSize = [1200, 800]
    render_view.Background = [0.32, 0.34, 0.43]  # Dark blue background

    # Display particles colored by extraction state
    particle_display = Show(particle_glyphs, render_view)
    particle_display.Representation = 'Surface'
    particle_display.ColorArrayName = ['POINTS', 'extraction_state']

    # Create custom color map for extraction (blue to brown)
    extraction_lut = GetColorTransferFunction('extraction_state')
    extraction_lut.RGBPoints = [
        0.0, 0.23, 0.30, 0.75,  # Blue (unextracted)
        0.5, 0.86, 0.64, 0.20,  # Orange (partial)
        1.0, 0.45, 0.27, 0.07   # Brown (fully extracted)
    ]
    extraction_lut.ColorSpace = 'HSV'

    # Add color bar
    color_bar = GetScalarBar(extraction_lut, render_view)
    color_bar.Title = 'Extraction State'
    color_bar.ComponentTitle = ''
    color_bar.Position = [0.85, 0.2]
    color_bar.ScalarBarLength = 0.6

    # Set up camera for good viewing angle
    camera = render_view.GetActiveCamera()
    camera.SetPosition(0.1, 0.1, 0.05)
    camera.SetFocalPoint(0.029, 0.015, 0.0)
    camera.SetViewUp(0.0, 0.0, 1.0)

    # Reset camera to fit data
    render_view.ResetCamera()

    # Add text annotations
    text_source = Text()
    text_source.Text = 'SproSim Espresso Extraction'
    text_display = Show(text_source, render_view)
    text_display.Color = [1.0, 1.0, 1.0]
    text_display.FontSize = 18
    text_display.Position = [0.02, 0.95]

    # Time annotation
    time_source = AnnotateTimeFilter(Input=reader)
    time_source.Format = 'Time: %.2f s'
    time_display = Show(time_source, render_view)
    time_display.Color = [1.0, 1.0, 1.0]
    time_display.FontSize = 14
    time_display.Position = [0.02, 0.05]

    print("Visualization setup complete!")

    # Set up animation
    animation_scene = GetAnimationScene()
    animation_scene.PlayMode = 'Snap To TimeSteps'

    if len(reader.TimestepValues) > 0:
        animation_scene.StartTime = reader.TimestepValues[0]
        animation_scene.EndTime = reader.TimestepValues[-1]
        print(f"Animation time range: {animation_scene.StartTime:.2f} - {animation_scene.EndTime:.2f} s")

    # Render the initial view
    Render()

    # Save screenshots if output directory specified
    if output_dir:
        save_visualization_outputs(render_view, animation_scene, output_dir)

    return render_view, particle_glyphs, reader

def save_visualization_outputs(render_view, animation_scene, output_dir):
    """Save screenshots and animation from the visualization."""

    os.makedirs(output_dir, exist_ok=True)

    # Save initial screenshot
    initial_screenshot = os.path.join(output_dir, 'espresso_initial.png')
    SaveScreenshot(initial_screenshot, render_view, ImageResolution=[1200, 800])
    print(f"Saved initial screenshot: {initial_screenshot}")

    # Save final frame screenshot
    animation_scene.GoToLast()
    Render()
    final_screenshot = os.path.join(output_dir, 'espresso_final.png')
    SaveScreenshot(final_screenshot, render_view, ImageResolution=[1200, 800])
    print(f"Saved final screenshot: {final_screenshot}")

    # Save animation
    animation_file = os.path.join(output_dir, 'espresso_animation.avi')
    SaveAnimation(animation_file, render_view,
                  ImageResolution=[1200, 800],
                  FrameRate=10)
    print(f"Saved animation: {animation_file}")

    # Reset to beginning
    animation_scene.GoToFirst()
    Render()

def create_flow_visualization(reader, render_view):
    """Add flow field visualization to existing setup."""

    print("Adding flow field visualization...")

    # Load flow field data (assuming separate flow files exist)
    # This would need to be adapted based on your actual flow data structure

    # For now, create a simple vector field representation
    # This is a placeholder - you'd need to load your actual flow VTS files

    # Create arrows for velocity vectors
    arrow_source = ArrowSource()
    arrow_source.TipRadius = 0.1
    arrow_source.ShaftRadius = 0.03

    # Note: This is simplified - you'd need to properly load and process
    # the flow field data from your VTS files

    print("Flow visualization setup complete!")

def main():
    """Main function for command line usage."""

    if len(sys.argv) < 2:
        print("Usage: pvpython visualize_espresso.py <path_to_brewing_simulation.pvd> [output_dir]")
        print("Example: pvpython visualize_espresso.py ./demo_output/brewing_simulation.pvd ./visualization_output")
        return

    pvd_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None

    if not os.path.exists(pvd_file):
        print(f"Error: File not found: {pvd_file}")
        return

    print("=== SproSim Espresso Visualization ===")
    print(f"Input file: {pvd_file}")
    if output_dir:
        print(f"Output directory: {output_dir}")

    try:
        render_view, glyphs, reader = setup_espresso_visualization(pvd_file, output_dir)

        print("\nVisualization complete!")
        print("Controls:")
        print("  - Use animation controls to play through time")
        print("  - Mouse: rotate (left), pan (middle), zoom (right)")
        print("  - Color bar shows extraction progress (blue → brown)")

        # Keep the visualization interactive
        if output_dir:
            print(f"\nOutputs saved to: {output_dir}")

    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()

# Interactive functions for ParaView Python shell
def quick_espresso_viz(pvd_file="brewing_simulation.pvd"):
    """Quick setup function for interactive use in ParaView."""
    return setup_espresso_visualization(pvd_file)

def export_images(output_dir="./paraview_output"):
    """Export current visualization as images."""
    render_view = GetActiveView()
    animation_scene = GetAnimationScene()
    if render_view:
        save_visualization_outputs(render_view, animation_scene, output_dir)
    else:
        print("No active view found. Run quick_espresso_viz() first.")

# Example usage instructions
USAGE_EXAMPLES = """
=== Usage Examples ===

1. Command line:
   pvpython visualize_espresso.py ./demo_output/brewing_simulation.pvd

2. In ParaView Python shell:
   exec(open('visualize_espresso.py').read())
   quick_espresso_viz('./demo_output/brewing_simulation.pvd')

3. Export images:
   export_images('./my_output_folder')

4. Custom visualization:
   render_view, glyphs, reader = setup_espresso_visualization('data.pvd')
   # Now customize as needed...
"""

if __name__ == '__main__':
    main()
else:
    # Print usage when imported
    print("SproSim ParaView visualization script loaded.")
    print("Run quick_espresso_viz('path/to/brewing_simulation.pvd') to start.")
