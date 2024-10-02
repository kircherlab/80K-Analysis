import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import io

# Create a simple plot
fig, ax = plt.subplots()
ax.plot([1, 2, 3], [1, 4, 9])

# Remove the default y-label
ax.set_ylabel("")

# Save the plot as an SVG in memory
f = io.StringIO()
plt.savefig(f, format="svg")
svg = f.getvalue()
f.close()

# Define the URL and the label you want to display
link = "https://www.example.com"
y_label = "Click here"

# Get the position where the y-label should be placed
trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transAxes)
y_label_x = -0.1  # This might need adjustment
y_label_y = 0.5  # Center of y-axis

# Create the SVG link element manually
svg_link = f'<text x="{y_label_x}" y="{y_label_y}" transform="scale(100,-100)" fill="black" ' \
           f'font-family="sans-serif" font-size="12">' \
           f'<a xlink:href="{link}" target="_blank">{y_label}</a>' \
           f'</text>'

# Insert the clickable y-label into the SVG content
insertion_point = svg.find('</svg>')
svg = svg[:insertion_point] + svg_link + svg[insertion_point:]

# Save the modified SVG to file
with open("plot_with_link.svg", "w") as f:
    f.write(svg)

# The resulting SVG will now have a clickable link as the y-label
