# Generates man/figures/camel-surface.png, the figure shown in README.md.
#
# Run from the package root with:  Rscript tools/make-camel-figure.R

library(OSCARS)

# ---- The objective: Branin's six-hump camel function ----
camel <- function(par) {
  x <- par[1]
  y <- par[2]
  4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4 - y^2)
}

# ---- Minimize over the box (-5,5) ----
set.seed(20260823)
out <- oscars(camel, n = 2, lwr = c(-5,-5), upr = c(5,5),
              controls = oscars.control(infol = 0))

# ---- The six minimizers that get marked ----
# Two global minimizers at (-+0.0898, +-0.7127), f = -1.0316, and four local
# ones: (-+1.7036, +-0.7961), f = -0.2155, and (+-1.6071, +-0.5687),
# f = 2.1043.  Seeds below are polished so the marks and the printed function
# values agree with the surface rather than with rounded literals.
refine <- function(p) {
  o <- optim(p, camel, method = "BFGS", control = list(reltol = 1e-15))
  list(par = o$par, value = o$value)
}
found  <- list(par = out$par, value = out$value)   # the one oscars returned
other  <- refine(-out$par)                         # the other global minimizer
locals <- lapply(list(c(-1.7036,  0.7961), c( 1.6071,  0.5687),
                      c( 1.7036, -0.7961), c(-1.6071, -0.5687)), refine)

localXY <- do.call(rbind, lapply(locals, `[[`, "par"))

# ---- The plotted window --------------------------------------------------
# Slightly wider than tall so all six minimizers sit clear of the edges.  The
# panel is sized below to keep one x unit the same length as one y unit.
xlim <- c(-2.25, 2.25)
ylim <- c(-2.00, 2.00)

x <- seq(xlim[1], xlim[2], length.out = 1351)
y <- seq(ylim[1], ylim[2], length.out = 1201)
z <- outer(x, y, function(x, y) {
  4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4 - y^2)
})

# The camel function spans roughly -1 to 60 on this window, so the surface is
# shaded on a log scale shifted to start at the global minimum.  On a linear
# scale the six humps flatten into one another.
fmin  <- found$value
tform <- function(f) log10(f - fmin + 1)
zt    <- tform(z)

pal  <- hcl.colors(512, "Viridis")
brks <- seq(min(zt), max(zt), length.out = length(pal) + 1)

# ---- Label helper --------------------------------------------------------
# A text box with a leader line running back to the point it describes.  The
# box is sized from the rendered text.  (lx, ly) is the centre of the box.
callout <- function(px, py, lx, ly, lines, cex = 0.8) {
  w    <- max(strwidth(lines, cex = cex, font = 2))
  h    <- strheight("Ag", cex = cex) * 1.75
  n    <- length(lines)
  padx <- strwidth("n", cex = cex)
  x0 <- lx - w/2 - padx
  x1 <- lx + w/2 + padx
  y0 <- ly - n*h/2 - 0.15*h
  y1 <- ly + n*h/2 + 0.15*h

  # leader runs to the nearest point on the box edge
  segments(px, py, min(max(px, x0), x1), min(max(py, y0), y1),
           col = "white", lwd = 1.4)
  rect(x0, y0, x1, y1, col = "#FFFFFFEE", border = "grey35", lwd = 1)
  for (i in seq_len(n)) {
    text(lx, ly + ((n - 1)/2 - (i - 1)) * h, lines[i], cex = cex,
         font = if (i == 1) 2 else 1,
         col  = if (i == 1) "black" else "grey25")
  }
}

# Labels sit along the top and bottom edges, in the same left/centre/right
# order as the points they describe, so every leader line runs straight down
# or straight up and none of them cross.
anchors <- list(found  = c( 0.00,  1.60),   # upper centre
                other  = c( 0.00, -1.62),   # lower centre
                local1 = c(-1.60,  1.60),   # upper left
                local2 = c( 1.60,  1.60),   # upper right
                local3 = c( 1.60, -1.62),   # lower right
                local4 = c(-1.60, -1.62))   # lower left

# ---- Draw ----------------------------------------------------------------
# Sizes are in inches so that the image panel is exactly 4.5 x 4.0 units of
# data at 1.4 in per unit, i.e. square data units in a rectangular panel.
pan <- c(6.30, 5.60)                     # image panel, width x height
mai <- c(0.70, 0.70, 0.45, 0.15)         # bottom, left, top, right
reg <- c(pan[1] + mai[2] + mai[4], pan[2] + mai[1] + mai[3])
barw <- 1.30                             # width of the colour bar column

png("man/figures/camel-surface.png",
    width = round((reg[1] + barw) * 200), height = round(reg[2] * 200),
    res = 200, type = "cairo", bg = "white")
on.exit(dev.off())

layout(matrix(c(1, 2), nrow = 1), widths = lcm(c(reg[1], barw) * 2.54))
par(mai = mai, cex.axis = 0.95)

image(x, y, zt, col = pal, breaks = brks, useRaster = TRUE,
      xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
axis(1); axis(2, las = 1); box()
mtext(expression(x[1]), side = 1, line = 2.4, cex = 1.0)
mtext(expression(x[2]), side = 2, line = 2.4, cex = 1.0)

levs <- c(-1, -0.5, 0, 0.5, 1, 2, 4, 8, 16, 32)
contour(x, y, z, levels = levs[levs > min(z) & levs < max(z)],
        drawlabels = FALSE, add = TRUE, col = "#FFFFFF55", lwd = 0.6)

# ---- labels first, then markers, so the leader lines run under the points --
callout(found$par[1], found$par[2], anchors$found[1], anchors$found[2],
        c("OSCARS minimum",
          sprintf("(%.4f, %.4f)", found$par[1], found$par[2]),
          sprintf("f = %.4f", found$value)))
callout(other$par[1], other$par[2], anchors$other[1], anchors$other[2],
        sprintf("global (f = %.4f)", other$value))
for (i in seq_along(locals)) {
  callout(locals[[i]]$par[1], locals[[i]]$par[2],
          anchors[[paste0("local", i)]][1], anchors[[paste0("local", i)]][2],
          sprintf("local (f = %.4f)", locals[[i]]$value))
}

# filled = the minimizer oscars returned, open red = the other global one,
# open grey = the four local ones
points(localXY[, 1], localXY[, 2], pch = 21, cex = 1.35, lwd = 1.8,
       col = "grey25", bg = "white")
points(other$par[1], other$par[2], pch = 21, cex = 1.7, lwd = 2.4,
       col = "#E8112D", bg = "white")
points(found$par[1], found$par[2], pch = 21, cex = 2.0, lwd = 2.4,
       col = "white", bg = "#E8112D")

# ---- Colour bar, labelled on the original f scale ----
par(mai = c(mai[1], 0.06, mai[3], 0.62))
plot(NA, xlim = c(0, 1), ylim = range(zt), xaxs = "i", yaxs = "i",
     axes = FALSE, xlab = "", ylab = "")
rasterImage(as.raster(matrix(rev(pal), ncol = 1)),
            0, min(zt), 1, max(zt), interpolate = TRUE)
ticks <- c(-1, 0, 1, 3, 10, 30, 60)
ticks <- ticks[tform(ticks) >= min(zt) & tform(ticks) <= max(zt)]
axis(4, at = tform(ticks), labels = ticks, las = 1, cex.axis = 0.8)
mtext("f(x)", side = 3, line = 1.0, cex = 0.85)
mtext("log scale", side = 3, line = 0.2, cex = 0.65, col = "grey35")
box()

cat("wrote man/figures/camel-surface.png\n")
cat(sprintf("oscars minimum %.6f at (%.6f, %.6f) in %d evaluations\n",
            out$value, out$par[1], out$par[2], out$evaluations))
