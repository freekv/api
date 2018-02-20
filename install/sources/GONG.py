import sunpy.map
class GONGMap(sunpy.map.GenericMap):
    def __init__(self, data, header, **kwargs):
        super(GONGMap, self).__init__(data, header, **kwargs)
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return header.get('telescop', '').startswith('NSO-GONG')
    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image."""
        return self.meta.get('wavelnth', 0)
