function save_figure(fig, outdir, basename)
% Save figure in multiple formats
    pngfile = fullfile(outdir, [basename, '.png']);
    pdffile = fullfile(outdir, [basename, '.pdf']);
    figfile = fullfile(outdir, [basename, '.fig']);

    try
        exportgraphics(fig, pngfile, 'Resolution', 300);
        exportgraphics(fig, pdffile, 'ContentType', 'vector');
    catch
        print(fig, pngfile, '-dpng', '-r300');
        print(fig, pdffile, '-dpdf', '-vector');
    end

    savefig(fig, figfile);
end